#ifndef ALILHCTAG_H
#define ALILHCTAG_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//-------------------------------------------------------------------------
//                          Class AliLHCTag
//   This is the class to deal with the tags for the LHC level
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <stdlib.h>
#include <Riostream.h>

#include "TObject.h"
#include "TClonesArray.h"

//______________________________________________________________________________
class AliLHCTag : public TObject
{
 private:
  Char_t   fLHCState[50];                 //LHC run conditions - comments
  Float_t  fLHCLuminosity;                //the value of the luminosity
  
 public:
  AliLHCTag();
  virtual ~AliLHCTag();
  
  void          SetLHCState(char *type) {strcpy(fLHCState,type);}
  void          SetLuminosity(Float_t lumin) {fLHCLuminosity = lumin;}
  void          SetLHCTag(Float_t lumin, char *type) {fLHCLuminosity = lumin; strcpy(fLHCState,type); }
  
  char         *GetLHCState() {return fLHCState;}
  Float_t       GetLuminosity() {return fLHCLuminosity;}
  
  ClassDef(AliLHCTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
