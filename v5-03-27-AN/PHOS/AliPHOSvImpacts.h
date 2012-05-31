#ifndef ALIPHOSVIMPACTS_H
#define ALIPHOSVIMPACTS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/* History of cvs commits:
 *
 * $Log$
 * Revision 1.9  2006/09/13 07:31:01  kharlov
 * Effective C++ corrections (T.Pocheptsov)
 *
 * Revision 1.8  2005/05/28 14:19:05  schutz
 * Compilation warnings fixed by T.P.
 *
 */

//_________________________________________________________________________
// Implementation version vImpacts of PHOS Manager class.
// This class inherits from v1 and adds impacts storing.
// Impacts stands for exact values of track coming to the detectors
// EMC, CPV or PPSD.
// Impacts are written to the same tree as hits are
// but in separate branches.
//*--                  
//*-- Author: Yuri Kharlov (IHEP, Protvino/SUBATECH, Nantes)

// --- ROOT system ---
class TList ;
class TLorentzVector ;


// --- AliRoot header files ---
#include "AliPHOSv1.h"

class AliPHOSvImpacts : public AliPHOSv1 {

public:

  AliPHOSvImpacts(void) ;
  AliPHOSvImpacts(const char *name, const char *title="") ;
  virtual ~AliPHOSvImpacts(void) ;

  virtual void   AddImpact(const char* detector, Int_t shunt, Int_t primary, Int_t track,
			   Int_t module, Int_t pid, TLorentzVector p, Float_t *xyz) ;
  virtual void   MakeBranch(Option_t *opt=" ");
  virtual void   ResetHits();
  virtual Int_t  IsVersion(void) const {
    // Gives the version number 
    return 1 ; 
  }
  virtual void   StepManager(void) ;                              
  virtual const TString Version(void)const{ 
    // returns the version number 
    return TString("vImpacts") ; 
  }

protected:

  TList * fEMCImpacts;         // Array of impacts in EMC modules
  TList * fCPVImpacts;         // Array of impacts in CPV (IHEP) modules
  TList * fPPSDImpacts;        // Array of impacts in PPSD modules

  Int_t  fNEMCImpacts[5];      // Number of EMC impacts per module
  Int_t  fNCPVImpacts[5];      // Number of CPV impacts per module
  Int_t  fNPPSDImpacts[5];     // Number of PPSD impacts per module

 private:
  AliPHOSvImpacts(AliPHOSvImpacts & phos);
  AliPHOSvImpacts & operator = (const AliPHOSvImpacts & /*rvalue*/);

  ClassDef(AliPHOSvImpacts,1)  // Implementation of PHOS manager class for layout EMC+PPSD

};

#endif // AliPHOSVIMPACTS_H
