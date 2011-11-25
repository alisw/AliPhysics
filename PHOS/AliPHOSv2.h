#ifndef ALIPHOSV2_H
#define ALIPHOSV2_H
/* Copyright(c) 1998-1999-2000, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.19  2005/07/01 20:01:36  kharlov
 * Warning fix on AddHit in gcc 3.4.2
 *
 * Revision 1.18  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Version of AliPHOSv0 which keeps all hits in TreeH
// I mean real hits not cumulated hits
//  This version is NOT recommended for Reconstruction analysis
//                  
//*-- Author: Gines MARTINEZ (SUBATECH)

// --- ROOT system ---

// --- AliRoot header files ---
#include "AliPHOSv1.h"

class AliPHOSv2 : public AliPHOSv1 {

public:

  AliPHOSv2(void) ;
  AliPHOSv2(const char *name, const char *title="") ;
  virtual ~AliPHOSv2(void) ;

  using AliPHOSv1::AddHit;
  virtual void    AddHit( Int_t shunt, Int_t primary, Int_t id, Float_t *hits); 
  virtual Int_t   IsVersion(void) const { 
    // Gives the version number 
    return 2 ; 
  }
  virtual const TString Version(void)const { 
    // returns the version number 
    return TString("v2") ; 
  }

private:

  AliPHOSv2(AliPHOSv2 & phos);
  AliPHOSv2 & operator = (const AliPHOSv2 & /*phos*/);

  ClassDef(AliPHOSv2,1)  // Class AliPHOSv0 which allows to write ond disk al the information of the hits. 

};

#endif // AliPHOSV2_H
