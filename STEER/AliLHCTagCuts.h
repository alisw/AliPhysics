#ifndef ALILHCTAGCUTS_H
#define ALILHCTAGCUTS_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliLHCTagCuts
//              This is the class for the cuts in run tags
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>

class AliLHCTag;

//___________________________________________________________________________
class AliLHCTagCuts : public TObject {
 public:
  AliLHCTagCuts();
  ~AliLHCTagCuts();
  void Reset();
  
 //____________________________________________________//
  void SetLHCState(TString state) {fLHCState = state; fLHCStateFlag = kTRUE;}
  void SetLHCLuminosityRange(Float_t low, Float_t high) {fLHCLuminosityMin = low; fLHCLuminosityMax = high; fLHCLuminosityFlag = kTRUE;}
 
  Bool_t IsAccepted(AliLHCTag *lhcTag) const;

  //____________________________________________________//
 private:
  TString fLHCState;              //LHC State
  Bool_t  fLHCStateFlag;          //Shows whether this cut is used or
  Float_t fLHCLuminosityMin;      //LHC luminosity - min
  Float_t fLHCLuminosityMax;      //LHC luminosity - max
  Bool_t  fLHCLuminosityFlag;     //Shows whether this cut is used or

  ClassDef(AliLHCTagCuts, 1)
};

#endif
