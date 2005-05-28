#ifndef ALIPHOSCPVRECPOINT_H
#define ALIPHOSCPVRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 */

//_________________________________________________________________________
//  RecPoint implementation for PHOS-CPV
//  An CpvRecPoint is a cluster of digits   
//*-- Author: Yuri Kharlov
//  (after Dmitri Peressounko (RRC KI & SUBATECH))
//  30 October 2000 
// --- ROOT system ---

//#include "TObject.h"
//#include "TArrayI.h"
 
// --- Standard library ---

// --- AliRoot header files ---

//#include "AliPHOSDigit.h"
#include "AliPHOSEmcRecPoint.h"
//#include "AliPHOSGeometry.h"

class AliPHOSCpvRecPoint : public AliPHOSEmcRecPoint  {

public:

  AliPHOSCpvRecPoint() ;
  AliPHOSCpvRecPoint(const char * opt) ;
  AliPHOSCpvRecPoint(const AliPHOSCpvRecPoint & rp) : AliPHOSEmcRecPoint(rp) {
    Fatal("cpy ctor", "not implemented") ; 
  } 
 
  virtual ~AliPHOSCpvRecPoint() ;  

  Int_t  Compare(const TObject * obj) const;                 // method for sorting  
  void   EvalAll(Float_t logWeight,TClonesArray * digits) ;
  void   EvalLocalPosition(Float_t logWeight,TClonesArray * digits ) ;  
  void   EvalClusterLengths(TClonesArray * digits) ;

  virtual void ExecuteEvent(Int_t event, Int_t px, Int_t py) /*const*/ ; 

  void   GetClusterLengths(Int_t &lengX, Int_t &lengZ) const {lengX = fLengX ;lengZ = fLengZ ;}
  Bool_t IsEmc(void) const {return kFALSE ; }        // tells that this is not a EMC
  Bool_t IsCPV(void) const {return kTRUE  ; }        // true if the recpoint is in CPV
  Bool_t IsSortable() const { return kTRUE ; }    // tells that this is a sortable object
  void   Print(const Option_t * = "") const ; 

  AliPHOSCpvRecPoint & operator = (const AliPHOSCpvRecPoint & /*rvalue*/)  {
   Fatal("operator =", "not implemented") ; return *this ; 
  }

 protected:

  Bool_t AreNeighbours(AliPHOSDigit * digit1, AliPHOSDigit * digit2 ) const ;

  Int_t    fLengX ;          // cluster length along x
  Int_t    fLengZ ;          // cluster length along z
  
  ClassDef(AliPHOSCpvRecPoint,1)  // CPV RecPoint (cluster)

};

#endif // AliPHOSCPVRECPOINT_H
