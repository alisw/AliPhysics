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

#include "TObject.h"
#include "TString.h"

//______________________________________________________________________________
class AliLHCTag : public TObject {
 public:
  AliLHCTag();
  virtual ~AliLHCTag();
  
  //____________________________________________________//
  void SetLHCState(TString type) {fLHCState = type;}
  void SetLuminosity(Float_t lumin) {fLHCLuminosity = lumin;}
  void SetLHCTag(Float_t lumin, TString type) {fLHCLuminosity = lumin; fLHCState = type; }
  
  //____________________________________________________//
  const char *GetLHCState() {return fLHCState.Data();}
  Float_t     GetLuminosity() const {return fLHCLuminosity;}
  
  //____________________________________________________//
 private:
  TString fLHCState;      //LHC run conditions - comments
  Float_t fLHCLuminosity; //the value of the luminosity
  
  ClassDef(AliLHCTag,1)  //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
