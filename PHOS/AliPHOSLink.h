#ifndef ALIPHOSLINK_H
#define ALIPHOSLINK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Algorithm class used only by AliPHOSTrackSegmentMaker 
//  Links recpoints
// into tracksegments                
//*-- Author: Dmitri Peressounko (SUBATECH)

// --- ROOT system ---

#include "TObject.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliPHOSLink : public  TObject{
  
public:
  
  AliPHOSLink( Float_t r, Int_t EMC, Int_t PPSD) ;  // ctor            
  virtual ~AliPHOSLink(){
    // dtor
  }
  Int_t   Compare(const TObject * obj) const;
  const Int_t   GetEmc(void) const { 
    // returns the index of EMC
    return fEmcN; 
  }
  const Int_t   GetPpsd(void) const { 
    // returns the index of PPSD
    return fPpsdN ; 
  }
  const Float_t GetR(void) const { 
    // returns the distance between EMC and PPSD
    return fR ; 
  } 
  Bool_t  IsSortable() const{ 
    // tells if this is a sortable object 
    return kTRUE ; 
  }
  
private:
  
  Int_t fEmcN ;  // Emc index
  Int_t fPpsdN ; // Ppsd index 
  Float_t fR ;   // Distance between EMC and PPSD RecPoints in a track segment 
  
  ClassDef(AliPHOSLink,1)  // Auxilliary algorithm class used by AliPHOSTrackSegmentMaker

};

#endif // AliPHOSLINK_H
