#ifndef ALIPHOSPPSDRECPOINT_H
#define ALIPHOSPPSDRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  A RecPoint (cluster) in the PPSD 
//  A PPSD RecPoint ends up to be a single digit
//  Oh yeah               
//*--  Yves Schutz (SUBATECH)

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSDigit.h"
#include "AliPHOSRecPoint.h"

class AliPHOSPpsdRecPoint : public AliPHOSRecPoint {

public:

  AliPHOSPpsdRecPoint() ;           // ctor   
  virtual ~AliPHOSPpsdRecPoint(){
    // dtor
  }
  virtual void AddDigit(AliPHOSDigit & digit, Float_t Energy) ;
  Int_t   Compare(const TObject * obj) const;                    // method to sort clusters

  void    EvalAll(Float_t LogWeight,TClonesArray * digits) ;
  void    EvalLocalPosition(Float_t LogWeight,TClonesArray * digits) ;   
                                                    // computes the position in the module of the cluster center
  void    EvalUp(TClonesArray * digits) ;              

  Int_t   GetMultiplicity(void)    const { return fMulDigit ;  } // returns the multiplicity of digits at 
                                                                 // the origin of this recpoint
  Int_t   GetMaximumMultiplicity() const { return   fMaxDigit ;} // returns the maximum allowed digit multiplicity 
  Bool_t  GetUp()      const {return fUp;     }    // true if cluster is in upper ppsd 
  Bool_t  IsEmc(void)  const {return kFALSE ; }    // tells that this is not a EMC
  Bool_t  IsSortable() const { return kTRUE ; }    // tells that this is a sortable object
  virtual void  Paint(Option_t * option="");
  void    Print(Option_t * option) ; 

private:

  Bool_t fUp ; 
  
  ClassDef(AliPHOSPpsdRecPoint,1)  // PPSD RecPoint

};

#endif // AliPHOSPPSDRECPOINT_H
