#ifndef ALIPHOSV0_HITS_H
#define ALIPHOSV0_HITS_H
/* Copyright(c) 1998-1999-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//_________________________________________________________________________
// Version of AliPHOSv0 which allows for keeping all hits in TreeH
//  This version is NOT recommended for Reconstruction analysis
//                  
//*-- Author: Gines MARTINEZ (SUBATECH)

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOSv0.h"
#include "AliPHOSReconstructioner.h"

class AliPHOSv0_hits : public AliPHOSv0 {

public:

  AliPHOSv0_hits(void) ;
  AliPHOSv0_hits(const char *name, const char *title="") ;
  virtual ~AliPHOSv0_hits(void) ;

  virtual void   AddHit( Int_t primary, Int_t id, Float_t *hits ) ; 
// adds a hit to the hit tree (any pre=digitalization is done here (so large root file !!) 
  void           FinishEvent(void) ;          // makes the digits from the hits 
protected:

  ClassDef(AliPHOSv0_hits,1)  // Class AliPHOSv0 which allows to write ond disk al the information of the hits. 

};

#endif // AliPHOSV0_HITS_H
