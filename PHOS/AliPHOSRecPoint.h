#ifndef ALIPHOSRECPOINT_H
#define ALIPHOSRECPOINT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_________________________________________________________________________
//  Base Class for PHOS Reconstructed Points  
//                  
//*-- Author: Gines Martinez (SUBATECH)

#include <assert.h>

// --- ROOT system ---

#include "TMarker.h"
#include "TGraph.h"
#include "TPaveText.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliRecPoint.h"
#include "AliPHOSDigit.h"
#include "AliPHOSGeometry.h"
class AliPHOSGeometry;

class AliPHOSRecPoint : public AliRecPoint {

 public:
  
  typedef TObjArray RecPointsList ; 

  AliPHOSRecPoint() ;                   // ctor         
  AliPHOSRecPoint(const AliPHOSRecPoint & rp) {
    // cpy ctor requested by Coding Convention 
    // but not yet needed
    assert(0==1) ; 
  } 
  
  virtual ~AliPHOSRecPoint(){
    // dtor
  }
  virtual  void   AddDigit(AliDigitNew &){
    // do not use this definition but the one below
    assert(0==1) ; 
  }
  virtual  void   AddDigit(AliPHOSDigit & digit, Float_t Energy) = 0 ; 
  virtual Int_t   Compare(TObject * obj) { 
    // check why this is done
    assert(0==1) ; return 1 ; 
  }   
  virtual Int_t   DistancetoPrimitive(Int_t px, Int_t py);
  virtual void    Draw(Option_t * option="") ;
  virtual void    ExecuteEvent(Int_t event, Int_t px, Int_t py) ;
  virtual Int_t   GetPHOSMod(void) ;
  virtual Int_t * GetPrimaries(Int_t & number) ;
  virtual Bool_t  IsEmc(void){
    // says that this is a EMC
    return kTRUE ;
  } 
  virtual Bool_t  IsSortable() const { 
    // tells that this is a sortable object
    return kTRUE ; 
  }  
  virtual void    Paint(Option_t * option="");
  virtual void    Print(Option_t * opt = "void") {
    // Print prototype
  } 

  AliPHOSRecPoint & operator = (const AliPHOSRecPoint & rvalue)  {
    // assignement operator requested by coding convention
    // but not needed
    assert(0==1) ;
    return *this ; 
  }

protected:
  
  Int_t fPHOSMod ;      // PHOS Module number in which the RecPoint is found
//    AliPHOSGeometry * fGeom ; // pointer to the PHOS geometry class
  
  ClassDef(AliPHOSRecPoint,1) // RecPoint for PHOS (Base Class)
 
};

#endif // AliPHOSRECPOINT_H
