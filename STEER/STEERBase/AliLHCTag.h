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
  AliLHCTag(const AliLHCTag &tag);
  virtual ~AliLHCTag();
  
  AliLHCTag &operator=(const AliLHCTag &tag);

  //____________________________________________________//
  //  void SetLHCTag(Float_t lumin, TString type) {fLHCLuminosity = lumin; fLHCState = type; }
  void UpdateFromRunTable(AliLHCTag &tag);

  void SetLHCState(TString type) {fLHCState = type;}
  void SetLuminosity(Float_t lumin) {fLHCLuminosity = lumin;}
  void SetNBunches(UShort_t nb) { fNBunches = nb; };
  void SetFillingScheme(TString sch) { fFillingScheme = sch; }
  void SetFillNo(Int_t fill) { fFillNo = fill; };
  void SetBeamEnergy(Float_t be) { fBeamEnergy = be; }
  void SetBunchIntensity(Float_t bi) { fBunchIntensity = bi; }
  
  
  //____________________________________________________//
  const char *GetLHCState() const {return fLHCState.Data();}
  Float_t     GetLuminosity() const {return fLHCLuminosity;}
  UShort_t    GetNBunches() const {return fNBunches; }
  TString     GetFillingScheme() const {return fFillingScheme; }
  Int_t       GetFillNo() const {return fFillNo; }
  Float_t     GetBeamEnergy() const {return fBeamEnergy; }
  Float_t     GetBunchIntensity() const {return fBunchIntensity; }

  //____________________________________________________//
 private:
  TString  fLHCState;      //LHC run conditions - comments
  Float_t  fLHCLuminosity; //the value of the luminosity
  UShort_t fNBunches;      //Number of bunches in beam
  TString  fFillingScheme; //Filling scheme name
  Int_t    fFillNo;        //Fill number
  Float_t  fBeamEnergy;    //Beam energy
  Float_t  fBunchIntensity;//Intensity per bunch

  ClassDef(AliLHCTag,2)   //(ClassName, ClassVersion)
};
//______________________________________________________________________________

#endif
