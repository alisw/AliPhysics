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
  Int_t   Compare(TObject * obj) ;                    // method to sort clusters

  Float_t GetDelta(void) {
    // returns the parameter used for sorting
    return fDelta ;
  }
  Int_t   GetMultiplicity(void) const { 
    // returns the multiplicity of digits at the origin of this recpoint
    return fMulDigit ; 
  } 
  Int_t   GetMaximumMultiplicity() { 
    // returns the maximum allowed digit multiplicity
    return   fMaxDigit ; 
  } 
  void    GetLocalPosition(TVector3 &LPos) ; // computes the position in the module of the cluster center 
  Float_t GetTotalEnergy(void) const { 
    // returns the amplitude for this recpoint 
    // in Ppsd EMC RecPoint Amp = Energy   
    return fAmp ; 
  }                             
  Bool_t  GetUp() ;               // true if cluster is in upper ppsd 
  Bool_t  IsEmc(void) {
    // tells that this is not a EMC
    return kFALSE ; 
  } 
  Bool_t  IsSortable() const { 
    // tells that this is a sortable object
    return kTRUE ; 
  }
  virtual void  Paint(Option_t * option="");
  void    Print(Option_t * opt = "void") ; 

private:

  Float_t        fDelta ;         // parameter used for sorting
  
  ClassDef(AliPHOSPpsdRecPoint,1)  // PPSD RecPoint

};

#endif // AliPHOSPPSDRECPOINT_H
