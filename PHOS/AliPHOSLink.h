#ifndef ALIPHOSLINK_H
#define ALIPHOSLINK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Short description                         //
//  Version SUBATECH                          //
//  Author Dmitri Peressounko   SUBATECH      //
//      comment: auxiliary class used   ONLY  //  
//               by AliPHOSTrackSegmentMaker  //
////////////////////////////////////////////////

// --- ROOT system ---

#include "TObject.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSLink : public  TObject{
  
public:
  
  AliPHOSLink( Float_t r, Int_t EMC, Int_t PPSD) ;  // ctor            
  virtual ~AliPHOSLink(){} // dtor
  
  Int_t   Compare(TObject * obj) ;
  Int_t   GetEmc(void) { return fEmcN; }
  Int_t   GetPpsd(void) { return fPpsdN ; }
  Float_t GetR(void) { return fR ; } 
  Bool_t  IsSortable() const{ return kTRUE ; }
  
private:
  
  Int_t fEmcN ;  // Emc index
  Int_t fPpsdN ; // Ppsd index 
  Float_t fR ;   // Distance 

public: 
  
  ClassDef(AliPHOSLink,1)  // description , version 1

};

#endif // AliPHOSLINK_H
