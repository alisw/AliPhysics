#ifndef ALIPHOSV2_H
#define ALIPHOSV2_H
/* Copyright(c) 1998-1999-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
//_________________________________________________________________________
// Version of AliPHOSv0 which keeps all hits in TreeH
// I mean real hits not cumulated hits
//  This version is NOT recommended for Reconstruction analysis
//                  
//*-- Author: Gines MARTINEZ (SUBATECH)

// --- ROOT system ---
#include "TClonesArray.h"

// --- AliRoot header files ---
#include "AliPHOSv1.h"
#include "AliPHOSReconstructioner.h"

class AliPHOSv2 : public AliPHOSv1 {

public:

  AliPHOSv2(void) ;
  AliPHOSv2(const char *name, const char *title="") ;
  virtual ~AliPHOSv2(void) ;

  virtual void    AddHit( Int_t shunt, Int_t primary, Int_t track, Int_t id, Float_t *hits, Int_t pid ) ; 
  virtual Int_t   IsVersion(void) const { 
    // Gives the version number 
    return 2 ; 
  }
  virtual TString Version(void){ 
    // returns the version number 
    return TString("v2") ; 
  }

protected:

  ClassDef(AliPHOSv2,1)  // Class AliPHOSv0 which allows to write ond disk al the information of the hits. 

};

#endif // AliPHOSV2_H
