#ifndef ALIDEBUGVOLUME_H
#define ALIDEBUGVOLUME_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-----------------------------------------------------------------------
//    Class to debug entry and exit from a volume
//    Used by AliLego class
//    Author: A.Morsch
//-----------------------------------------------------------------------

#include "TNamed.h"
class AliDebugVolume : public TNamed  {

public:
  AliDebugVolume();
  AliDebugVolume(const char *name, Int_t copy,
		 Float_t step, Float_t x, Float_t y, Float_t z, Int_t status);
  virtual ~AliDebugVolume(){}
  
  Int_t   CopyNumber() const {return fCopy;}
  Float_t Step()       const {return fStep;}
  Float_t X()          const {return fX;}  
  Float_t Y()          const {return fY;}
  Float_t Z()          const {return fZ;}
  const char*   Status()     const;
  
  
  Bool_t  IsVEqual(const char* name, Int_t copy) const;
private:
   Int_t      fCopy;             //!Volume copy number
   Float_t    fStep;             //!Stepsize to volume boundary
   Float_t    fX;                // x
   Float_t    fY;                // y
   Float_t    fZ;                // z of boundary crossing
   Int_t      fStatus;           // tracking status
   
  ClassDef(AliDebugVolume,1)      //Utility class to store volume information
                                  //during debugging 

};


#endif
